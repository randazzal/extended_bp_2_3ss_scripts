library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(ggeasy)

setwd("/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/")
tpms <- read.delim("all_tpms.txt", sep = "\t")
psi <- read.delim("megacombo_sizelimited.txt", sep = "\t")
psi1 <- psi[,1]
psi1 <- as.data.frame(psi1)
psi2 <- psi[,6:15]
psi <- cbind(psi1, psi2)
heatmap <- psi %>%
  pivot_longer(
    cols = c("IncLevel_mean_HH11", "IncLevel_mean_HH15", "IncLevel_mean_day5",
             "IncLevel_mean_day7", "IncLevel_mean_day9", "IncLevel_mean_HH11h",
             "IncLevel_mean_HH15h", "IncLevel_mean_day5h", "IncLevel_mean_day7h",
             "IncLevel_mean_day9h"),
    names_to = "name",
    values_to = "PSI"
  ) %>%
  mutate(
    Tissue = if_else(grepl("h$", name), "H", "B"),
    Stage  = gsub("\\D", "", gsub("^IncLevel_mean_", "", name))  # keep digits only
  )

gene_id <- read.delim("/data2/lackey_lab/DownloadedSequenceData/randazza/new_chick/gene_id_name.txt", sep = "\t", header = FALSE)

tpms <- merge(tpms, gene_id, by.x = "IDs", by.y = "V1")

limited_long <- tpms %>%
  pivot_longer(cols = c("Male_5B1", "Female_5B1", "Female_5B2", "Male_5B2", "Female_5B3", "Male_5B3",
                        "Male_5H1", "Female_5H1", "Female_5H2", "Female_5H3", "Male_5H2", "Male_5H3",
                        "Female_7B1", "Male_7B1", "Female_7B2", "Female_7B3", "Male_7B2", "Male_7B3",
                        "Female_7H1", "Female_7H2", "Male_7H1", "Female_7H3", "Male_7H2", "Male_7H3",
                        "Female_9B1", "Male_9B1", "Male_9B2", "Female_9B2", "Male_9B3", "Female_9B3",
                        "Female_9H1", "Female_9H2", "Female_9H3", "Male_9H1", "Male_9H2", "Female_9H4", "Male_9H3",
                        "Male_11B1", "Male_11B2", "Female_11B1", "Female_11B2", "Male_11B3", "Male_11B4", "Female_11B3",
                        "Male_11H1", "Male_11H2", "Female_11H1", "Male_11H3", "Male_11H4", "Male_11H5", "Female_11H2", "Female_11H3",
                        "Female_15B1", "Female_15B2", "Male_15B1", "Male_15B2", "Male_15B3", "Female_15B3", "Male_15B4",
                        "Female_15H1", "Female_15H2", "Female_15H3", "Male_15H1", "Male_15H2", "Female_15H4", "Male_15H3"),
               names_to = "Sample", values_to = "TPMs")
tpm_long <- limited_long %>%
  mutate(
    Sex = str_extract(Sample, "Male|Female"),
    Stage = str_extract(Sample, "(?<=_)[0-9]+"),
    Tissue = sub(".*_[0-9]+([A-Z])[0-9]+$", "\\1", Sample),
    Replicate = str_extract(Sample, "[0-9]+$")
  )
tpm_summary <- tpm_long %>%
  group_by(IDs, V2, Tissue, Stage) %>%
  summarise(mean_TPM = mean(TPMs, na.rm = TRUE), .groups = "drop")
colnames(tpm_summary) <- c("IDs", "psi1", "Tissue", "Stage", "mean_TPM")
combo <- merge(tpm_summary, heatmap, by = c("Tissue", "Stage", "psi1"))

corr_val <- cor(combo$PSI, combo$mean_TPM, 
                use = "complete.obs", method = "pearson")
combo %>%
  ggplot(aes(x = PSI, y = mean_TPM)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "PSI vs TPM",
    x = "PSI (%)",
    y = "Mean TPM",
    subtitle = paste("Pearson Correlation =", round(corr_val, 3))) +
  ggeasy::easy_center_title() +
  geom_smooth(method = "lm", se = FALSE, color = "orange", linetype = "dashed")
