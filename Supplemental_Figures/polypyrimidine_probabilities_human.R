library(dplyr)
library(tidyr)
library(stringr)
library(Biostrings)
library(readxl)
library(ggplot2)
library(fs)
library(purrr)

setwd("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/")
################################# A3'ss #####################################
set.seed(123)
neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_neuro_A.fasta")
# Function to calculate long-format probability matrix
calc_pfm_long <- function(seqs) {
  pfm_counts <- consensusMatrix(seqs)
  pfm_counts <- pfm_counts[c("A", "C", "G", "T"), , drop = FALSE]
  n_seqs <- length(seqs)
  df <- as.data.frame(pfm_counts)
  df$nucleotide <- rownames(df)
  rownames(df) <- NULL
  df_long <- df %>%
    pivot_longer(cols = -nucleotide,
                 names_to = "position",
                 values_to = "count"
                 ) %>%
    mutate(position = as.numeric(sub("V", "", position)),
           total = n_seqs,
           probability = count / n_seqs
           )
  return(df_long)
}

micro_neuro <- calc_pfm_long(neuro_A) %>%
  mutate(group = "Neural Microexons")
summary_neuro <- calc_pfm_long(neuro_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_neuro$group <- "Neural Microexons"

nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_non_neuro_A.fasta")
micro_nn <- calc_pfm_long(nn_A) %>%
  mutate(group = "Non-Neural Microexons")
summary_nn <- calc_pfm_long(nn_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_nn$group <- "Non-neural Microexon"

avg_nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_non_neuro_A.fasta")
mid_nn <- calc_pfm_long(avg_nn_A) %>%
  mutate(group = "Non-Neural Midexons")
summary_avg_nn <- calc_pfm_long(avg_nn_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_avg_nn$group <- "Non-neural Midexon"

avg_neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_neuro_A.fasta")
mid_neuro <- calc_pfm_long(avg_neuro_A) %>%
  mutate(group = "Neural Midexons")
summary_avg_neuro <- calc_pfm_long(avg_neuro_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_avg_neuro$group <- "Neural Midexon"

mega_summary <- rbind(summary_neuro, summary_nn, summary_avg_neuro, summary_avg_nn)
mega_summary %>%
  filter(position > 399,
         position < 453) %>% 
  ggplot(aes(x = position, y = mean_prob, color = nucleotide)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Bootstrapped PFM (n = 100 each)", y = "Probability") +
  facet_grid(rows = vars(group))

## U C combo
avg_neu_polyU <- summary_avg_neuro %>%
  filter(nucleotide %in% c("T"))
avg_neu_polyU <- avg_neu_polyU[ ,c("position", "mean_prob")]
colnames(avg_neu_polyU) <- c("position", "mean_prob_T")
avg_neu_polyC <- summary_avg_neuro %>%
  filter(nucleotide %in% c("C"))
avg_neu_polyC <- avg_neu_polyC[ ,c("position", "mean_prob")]
colnames(avg_neu_polyC) <- c("position", "mean_prob_C")
avg_neu_UC <- merge(avg_neu_polyU, avg_neu_polyC, by = "position")
avg_neu_UC$UC <- avg_neu_UC$mean_prob_T + avg_neu_UC$mean_prob_C
avg_neu_UC$group <- "Neural Midexons"

avg_nn_polyU <- summary_avg_nn %>%
  filter(nucleotide %in% c("T"))
avg_nn_polyU <- avg_nn_polyU[ ,c("position", "mean_prob")]
colnames(avg_nn_polyU) <- c("position", "mean_prob_T")
avg_nn_polyC <- summary_avg_nn %>%
  filter(nucleotide %in% c("C"))
avg_nn_polyC <- avg_nn_polyC[ ,c("position", "mean_prob")]
colnames(avg_nn_polyC) <- c("position", "mean_prob_C")
avg_nn_UC <- merge(avg_nn_polyU, avg_nn_polyC, by = "position")
avg_nn_UC$UC <- avg_nn_UC$mean_prob_T + avg_nn_UC$mean_prob_C
avg_nn_UC$group <- "Non-Neural Midexons"

mic_neu_polyU <- summary_neuro %>%
  filter(nucleotide %in% c("T"))
mic_neu_polyU <- mic_neu_polyU[ ,c("position", "mean_prob")]
colnames(mic_neu_polyU) <- c("position", "mean_prob_T")
mic_neu_polyC <- summary_neuro %>%
  filter(nucleotide %in% c("C"))
mic_neu_polyC <- mic_neu_polyC[ ,c("position", "mean_prob")]
colnames(mic_neu_polyC) <- c("position", "mean_prob_C")
mic_neu_UC <- merge(mic_neu_polyU, mic_neu_polyC, by = "position")
mic_neu_UC$UC <- mic_neu_UC$mean_prob_T + mic_neu_UC$mean_prob_C
mic_neu_UC$group <- "Neural Microexons"

mic_nn_polyU <- summary_nn %>%
  filter(nucleotide %in% c("T"))
mic_nn_polyU <- mic_nn_polyU[ ,c("position", "mean_prob")]
colnames(mic_nn_polyU) <- c("position", "mean_prob_T")
mic_nn_polyC <- summary_nn %>%
  filter(nucleotide %in% c("C"))
mic_nn_polyC <- mic_nn_polyC[ ,c("position", "mean_prob")]
colnames(mic_nn_polyC) <- c("position", "mean_prob_C")
mic_nn_UC <- merge(mic_nn_polyU, mic_nn_polyC, by = "position")
mic_nn_UC$UC <- mic_nn_UC$mean_prob_T + mic_nn_UC$mean_prob_C
mic_nn_UC$group <- "Non-Neural Microexon"

polypy <- rbind(mic_neu_UC, mic_nn_UC, avg_neu_UC, avg_nn_UC)
polypy %>%
  filter(position > 399,
         position < 453) %>% 
  ggplot(aes(x = position, y = UC, color = group)) +
  geom_line() +
  theme_minimal() +
  labs(title = "PolyPyrimidine Probability", y = "Probability") +
  scale_color_manual(values = c("#fdb863", "#5e3c99", "#e66101", "#b2abd2")) +
  ylim(0,1)


calc_pfm_long_counts <- function(seqs) {
  pfm_counts <- consensusMatrix(seqs)
  pfm_counts <- pfm_counts[c("A", "C", "G", "T"), , drop = FALSE]
  n_seqs <- length(seqs)
  
  df <- as.data.frame(pfm_counts)
  df$nucleotide <- rownames(df)
  rownames(df) <- NULL
  
  df_long <- df %>%
    pivot_longer(cols = -nucleotide,
                 names_to = "position",
                 values_to = "count") %>%
    mutate(position = as.numeric(sub("V", "", position)),
           pyrimidine = nucleotide %in% c("C", "T"))
  
  return(df_long)
}
neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_neuro_A.fasta")
micro_neuro <- calc_pfm_long_counts(neuro_A) %>%
  mutate(group = "Neural Microexons")
nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_non_neuro_A.fasta")
micro_nn <- calc_pfm_long_counts(nn_A) %>%
  mutate(group = "Non-Neural Microexons")
avg_nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_non_neuro_A.fasta")
mid_nn <- calc_pfm_long_counts(avg_nn_A) %>%
  mutate(group = "Non-Neural Midexons")
avg_neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_neuro_A.fasta")
mid_neuro <- calc_pfm_long_counts(avg_neuro_A) %>%
  mutate(group = "Neural Midexons")

polypy <- rbind(mid_neuro, mid_nn, micro_neuro, micro_nn)
compare_groups_corrected <- function(df, g1, g2) {
  df %>%
    filter(group %in% c(g1, g2)) %>%
    group_by(position) %>%
    summarise({
      table_counts <- cur_data() %>%
        group_by(group, pyrimidine) %>%
        summarise(total = sum(count), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = pyrimidine, values_from = total, values_fill = 0)
      
      mat <- as.matrix(table_counts[, c("TRUE", "FALSE")])
      rownames(mat) <- table_counts$group
      test <- fisher.test(mat)
      
      tibble(
        p_value = test$p.value,
        odds_ratio = test$estimate
      )
    }, .groups = "drop")
}

p_nmicro_navg <- compare_groups_corrected(polypy, "Neural Microexons", "Neural Midexons")
p_nmicro_nonm <- compare_groups_corrected(polypy, "Neural Microexons", "Non-Neural Microexons")
p_nmicro_nonavg <- compare_groups_corrected(polypy, "Neural Microexons", "Non-Neural Midexons")

# Merge results
p_values_df <- purrr::reduce(list(
  p_nmicro_navg, #p_value.x
  p_nmicro_nonm, #p_value.y
  p_nmicro_nonavg #p_value
), dplyr::full_join, by = "position")

# Adjust p-values for multiple comparisons 
p_values_df <- p_values_df %>%
  mutate(across(starts_with("p_value"), ~ p.adjust(.x, method = "fdr")))
p_values_df_fixed <- p_values_df %>%
  rename_with(~ paste0(., ".z"), .cols = c("p_value", "odds_ratio"))
Longer <- p_values_df_fixed %>%
  filter(position > 399, position < 453) %>%
  pivot_longer(
    cols = c(odds_ratio.x, odds_ratio.y, odds_ratio.z,
             p_value.x, p_value.y, p_value.z),
    names_to = c(".value", "comparison"),
    names_pattern = "(.*)\\.(.*)"
  ) %>%
  mutate(
    comparison = recode(comparison,
                        "x" = "Micro vs Avg",
                        "y" = "Micro vs NonM",
                        "z" = "Micro vs NonAvg"),
    log2OR = log2(odds_ratio),
    signif = p_value < 0.05
  )
ggplot(Longer, aes(x = position, y = log2OR, color = signif)) +
  geom_line(aes(group = 1), linewidth = 1) +           # trend line
  geom_point(size = 2) +                          # points for each position
  geom_hline(yintercept = 0, linetype = "dashed") +  # baseline = no enrichment
  scale_color_manual(values = c("grey70", "orange")) +  # non-significant vs significant
  facet_wrap(~comparison, ncol = 1) +            # one panel per comparison
  theme_minimal(base_size = 14) +
  labs(
    x = "Nucleotide position",
    y = "log2(Odds Ratio) (Pyrimidine Enrichment)",
    color = "FDR < 0.05",
    title = "Position-wise Pyrimidine Enrichment in Microexons"
  )

################################# A5'ss #######################################
neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_neuro_A5.fasta")
# Function to calculate long-format probability matrix
calc_pfm_long <- function(seqs) {
  pfm_counts <- consensusMatrix(seqs)
  pfm_counts <- pfm_counts[c("A", "C", "G", "T"), , drop = FALSE]
  n_seqs <- length(seqs)
  df <- as.data.frame(pfm_counts)
  df$nucleotide <- rownames(df)
  rownames(df) <- NULL
  df_long <- df %>%
    pivot_longer(cols = -nucleotide,
                 names_to = "position",
                 values_to = "count"
    ) %>%
    mutate(position = as.numeric(sub("V", "", position)),
           total = n_seqs,
           probability = count / n_seqs
    )
  return(df_long)
}

micro_neuro <- calc_pfm_long(neuro_A) %>%
  mutate(group = "Neural Microexons")
summary_neuro <- calc_pfm_long(neuro_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_neuro$group <- "Neural Microexons"

nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_non_neuro_A5.fasta")
micro_nn <- calc_pfm_long(nn_A) %>%
  mutate(group = "Non-Neural Microexons")
summary_nn <- calc_pfm_long(nn_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_nn$group <- "Non-neural Microexon"

avg_nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_non_neuro_A5.fasta")
mid_nn <- calc_pfm_long(avg_nn_A) %>%
  mutate(group = "Non-Neural Midexons")
summary_avg_nn <- calc_pfm_long(avg_nn_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_avg_nn$group <- "Non-neural Midexon"

avg_neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_neuro_A5.fasta")
mid_neuro <- calc_pfm_long(avg_neuro_A) %>%
  mutate(group = "Neural Midexons")
summary_avg_neuro <- calc_pfm_long(avg_neuro_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_avg_neuro$group <- "Neural Midexon"

mega_summary <- rbind(summary_neuro, summary_nn, summary_avg_neuro, summary_avg_nn)
mega_summary %>%
  filter(position >= 451,
         position < 502) %>% 
  ggplot(aes(x = position, y = mean_prob, color = nucleotide)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Bootstrapped PFM (n = 100 each)", y = "Probability") +
  facet_grid(rows = vars(group))

## U C combo
avg_neu_polyU <- summary_avg_neuro %>%
  filter(nucleotide %in% c("T"))
avg_neu_polyU <- avg_neu_polyU[ ,c("position", "mean_prob")]
colnames(avg_neu_polyU) <- c("position", "mean_prob_T")
avg_neu_polyC <- summary_avg_neuro %>%
  filter(nucleotide %in% c("C"))
avg_neu_polyC <- avg_neu_polyC[ ,c("position", "mean_prob")]
colnames(avg_neu_polyC) <- c("position", "mean_prob_C")
avg_neu_UC <- merge(avg_neu_polyU, avg_neu_polyC, by = "position")
avg_neu_UC$UC <- avg_neu_UC$mean_prob_T + avg_neu_UC$mean_prob_C
avg_neu_UC$group <- "Neural Midexons"

avg_nn_polyU <- summary_avg_nn %>%
  filter(nucleotide %in% c("T"))
avg_nn_polyU <- avg_nn_polyU[ ,c("position", "mean_prob")]
colnames(avg_nn_polyU) <- c("position", "mean_prob_T")
avg_nn_polyC <- summary_avg_nn %>%
  filter(nucleotide %in% c("C"))
avg_nn_polyC <- avg_nn_polyC[ ,c("position", "mean_prob")]
colnames(avg_nn_polyC) <- c("position", "mean_prob_C")
avg_nn_UC <- merge(avg_nn_polyU, avg_nn_polyC, by = "position")
avg_nn_UC$UC <- avg_nn_UC$mean_prob_T + avg_nn_UC$mean_prob_C
avg_nn_UC$group <- "Non-Neural Midexons"

mic_neu_polyU <- summary_neuro %>%
  filter(nucleotide %in% c("T"))
mic_neu_polyU <- mic_neu_polyU[ ,c("position", "mean_prob")]
colnames(mic_neu_polyU) <- c("position", "mean_prob_T")
mic_neu_polyC <- summary_neuro %>%
  filter(nucleotide %in% c("C"))
mic_neu_polyC <- mic_neu_polyC[ ,c("position", "mean_prob")]
colnames(mic_neu_polyC) <- c("position", "mean_prob_C")
mic_neu_UC <- merge(mic_neu_polyU, mic_neu_polyC, by = "position")
mic_neu_UC$UC <- mic_neu_UC$mean_prob_T + mic_neu_UC$mean_prob_C
mic_neu_UC$group <- "Neural Microexons"

mic_nn_polyU <- summary_nn %>%
  filter(nucleotide %in% c("T"))
mic_nn_polyU <- mic_nn_polyU[ ,c("position", "mean_prob")]
colnames(mic_nn_polyU) <- c("position", "mean_prob_T")
mic_nn_polyC <- summary_nn %>%
  filter(nucleotide %in% c("C"))
mic_nn_polyC <- mic_nn_polyC[ ,c("position", "mean_prob")]
colnames(mic_nn_polyC) <- c("position", "mean_prob_C")
mic_nn_UC <- merge(mic_nn_polyU, mic_nn_polyC, by = "position")
mic_nn_UC$UC <- mic_nn_UC$mean_prob_T + mic_nn_UC$mean_prob_C
mic_nn_UC$group <- "Non-Neural Microexon"

polypy <- rbind(mic_neu_UC, mic_nn_UC, avg_neu_UC, avg_nn_UC)
polypy %>%
  filter(position >= 451,
         position < 502) %>% 
  ggplot(aes(x = position, y = UC, color = group)) +
  geom_line() +
  theme_minimal() +
  labs(title = "PolyPyrimidine Probability", y = "Probability") +
  scale_color_manual(values = c("#fdb863", "#5e3c99", "#e66101", "#b2abd2")) +
  ylim(0,1)

neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_neuro_A5.fasta")
micro_neuro <- calc_pfm_long_counts(neuro_A) %>%
  mutate(group = "Neural Microexons")
nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_non_neuro_A5.fasta")
micro_nn <- calc_pfm_long_counts(nn_A) %>%
  mutate(group = "Non-Neural Microexons")
avg_nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_non_neuro_A5.fasta")
mid_nn <- calc_pfm_long_counts(avg_nn_A) %>%
  mutate(group = "Non-Neural Midexons")
avg_neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_neuro_A5.fasta")
mid_neuro <- calc_pfm_long_counts(avg_neuro_A) %>%
  mutate(group = "Neural Midexons")

polypy <- rbind(mid_neuro, mid_nn, micro_neuro, micro_nn)
compare_groups_corrected <- function(df, g1, g2) {
  df %>%
    filter(group %in% c(g1, g2)) %>%
    group_by(position) %>%
    summarise({
      table_counts <- cur_data() %>%
        group_by(group, pyrimidine) %>%
        summarise(total = sum(count), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = pyrimidine, values_from = total, values_fill = 0)
      
      mat <- as.matrix(table_counts[, c("TRUE", "FALSE")])
      rownames(mat) <- table_counts$group
      test <- fisher.test(mat)
      
      tibble(
        p_value = test$p.value,
        odds_ratio = test$estimate
      )
    }, .groups = "drop")
}

p_nmicro_navg <- compare_groups_corrected(polypy, "Neural Microexons", "Neural Midexons")
p_nmicro_nonm <- compare_groups_corrected(polypy, "Neural Microexons", "Non-Neural Microexons")
p_nmicro_nonavg <- compare_groups_corrected(polypy, "Neural Microexons", "Non-Neural Midexons")

# Merge results
p_values_df <- purrr::reduce(list(
  p_nmicro_navg, #p_value.x
  p_nmicro_nonm, #p_value.y
  p_nmicro_nonavg #p_value
), dplyr::full_join, by = "position")

# Adjust p-values for multiple comparisons 
p_values_df <- p_values_df %>%
  mutate(across(starts_with("p_value"), ~ p.adjust(.x, method = "fdr")))
p_values_df_fixed <- p_values_df %>%
  rename_with(~ paste0(., ".z"), .cols = c("p_value", "odds_ratio"))
Longer <- p_values_df_fixed %>%
  filter(position >= 451,
        position < 502) %>%
  pivot_longer(
    cols = c(odds_ratio.x, odds_ratio.y, odds_ratio.z,
             p_value.x, p_value.y, p_value.z),
    names_to = c(".value", "comparison"),
    names_pattern = "(.*)\\.(.*)"
  ) %>%
  mutate(
    comparison = recode(comparison,
                        "x" = "Micro vs Avg",
                        "y" = "Micro vs NonM",
                        "z" = "Micro vs NonAvg"),
    log2OR = log2(odds_ratio),
    signif = p_value < 0.05
  )
ggplot(Longer, aes(x = position, y = log2OR, color = signif)) +
  geom_line(aes(group = 1), linewidth = 1) +           # trend line
  geom_point(size = 2) +                          # points for each position
  geom_hline(yintercept = 0, linetype = "dashed") +  # baseline = no enrichment
  scale_color_manual(values = c("grey70", "orange")) +  # non-significant vs significant
  facet_wrap(~comparison, ncol = 1) +            # one panel per comparison
  theme_minimal(base_size = 14) +
  labs(
    x = "Nucleotide position",
    y = "log2(Odds Ratio) (Pyrimidine Enrichment)",
    color = "FDR < 0.05",
    title = "Position-wise Pyrimidine Enrichment in Microexons"
  )

################################# C23'ss #####################################
neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_neuro_C2.fasta")
# Function to calculate long-format probability matrix
calc_pfm_long <- function(seqs) {
  pfm_counts <- consensusMatrix(seqs)
  pfm_counts <- pfm_counts[c("A", "C", "G", "T"), , drop = FALSE]
  n_seqs <- length(seqs)
  df <- as.data.frame(pfm_counts)
  df$nucleotide <- rownames(df)
  rownames(df) <- NULL
  df_long <- df %>%
    pivot_longer(cols = -nucleotide,
                 names_to = "position",
                 values_to = "count"
    ) %>%
    mutate(position = as.numeric(sub("V", "", position)),
           total = n_seqs,
           probability = count / n_seqs
    )
  return(df_long)
}

micro_neuro <- calc_pfm_long(neuro_A) %>%
  mutate(group = "Neural Microexons")
summary_neuro <- calc_pfm_long(neuro_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_neuro$group <- "Neural Microexons"

nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_non_neuro_C2.fasta")
micro_nn <- calc_pfm_long(nn_A) %>%
  mutate(group = "Non-Neural Microexons")
summary_nn <- calc_pfm_long(nn_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_nn$group <- "Non-neural Microexon"

avg_nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_non_neuro_C2.fasta")
mid_nn <- calc_pfm_long(avg_nn_A) %>%
  mutate(group = "Non-Neural Midexons")
summary_avg_nn <- calc_pfm_long(avg_nn_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_avg_nn$group <- "Non-neural Midexon"

avg_neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_neuro_C2.fasta")
mid_neuro <- calc_pfm_long(avg_neuro_A) %>%
  mutate(group = "Neural Midexons")
summary_avg_neuro <- calc_pfm_long(avg_neuro_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_avg_neuro$group <- "Neural Midexon"

mega_summary <- rbind(summary_neuro, summary_nn, summary_avg_neuro, summary_avg_nn)
mega_summary %>%
  filter(position > 649,
         position < 703) %>% 
  ggplot(aes(x = position, y = mean_prob, color = nucleotide)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Bootstrapped PFM (n = 100 each)", y = "Probability") +
  facet_grid(rows = vars(group))

## U C combo
avg_neu_polyU <- summary_avg_neuro %>%
  filter(nucleotide %in% c("T"))
avg_neu_polyU <- avg_neu_polyU[ ,c("position", "mean_prob")]
colnames(avg_neu_polyU) <- c("position", "mean_prob_T")
avg_neu_polyC <- summary_avg_neuro %>%
  filter(nucleotide %in% c("C"))
avg_neu_polyC <- avg_neu_polyC[ ,c("position", "mean_prob")]
colnames(avg_neu_polyC) <- c("position", "mean_prob_C")
avg_neu_UC <- merge(avg_neu_polyU, avg_neu_polyC, by = "position")
avg_neu_UC$UC <- avg_neu_UC$mean_prob_T + avg_neu_UC$mean_prob_C
avg_neu_UC$group <- "Neural Midexons"

avg_nn_polyU <- summary_avg_nn %>%
  filter(nucleotide %in% c("T"))
avg_nn_polyU <- avg_nn_polyU[ ,c("position", "mean_prob")]
colnames(avg_nn_polyU) <- c("position", "mean_prob_T")
avg_nn_polyC <- summary_avg_nn %>%
  filter(nucleotide %in% c("C"))
avg_nn_polyC <- avg_nn_polyC[ ,c("position", "mean_prob")]
colnames(avg_nn_polyC) <- c("position", "mean_prob_C")
avg_nn_UC <- merge(avg_nn_polyU, avg_nn_polyC, by = "position")
avg_nn_UC$UC <- avg_nn_UC$mean_prob_T + avg_nn_UC$mean_prob_C
avg_nn_UC$group <- "Non-Neural Midexons"

mic_neu_polyU <- summary_neuro %>%
  filter(nucleotide %in% c("T"))
mic_neu_polyU <- mic_neu_polyU[ ,c("position", "mean_prob")]
colnames(mic_neu_polyU) <- c("position", "mean_prob_T")
mic_neu_polyC <- summary_neuro %>%
  filter(nucleotide %in% c("C"))
mic_neu_polyC <- mic_neu_polyC[ ,c("position", "mean_prob")]
colnames(mic_neu_polyC) <- c("position", "mean_prob_C")
mic_neu_UC <- merge(mic_neu_polyU, mic_neu_polyC, by = "position")
mic_neu_UC$UC <- mic_neu_UC$mean_prob_T + mic_neu_UC$mean_prob_C
mic_neu_UC$group <- "Neural Microexons"

mic_nn_polyU <- summary_nn %>%
  filter(nucleotide %in% c("T"))
mic_nn_polyU <- mic_nn_polyU[ ,c("position", "mean_prob")]
colnames(mic_nn_polyU) <- c("position", "mean_prob_T")
mic_nn_polyC <- summary_nn %>%
  filter(nucleotide %in% c("C"))
mic_nn_polyC <- mic_nn_polyC[ ,c("position", "mean_prob")]
colnames(mic_nn_polyC) <- c("position", "mean_prob_C")
mic_nn_UC <- merge(mic_nn_polyU, mic_nn_polyC, by = "position")
mic_nn_UC$UC <- mic_nn_UC$mean_prob_T + mic_nn_UC$mean_prob_C
mic_nn_UC$group <- "Non-Neural Microexon"

polypy <- rbind(mic_neu_UC, mic_nn_UC, avg_neu_UC, avg_nn_UC)
polypy %>%
  filter(position > 649,
         position < 703) %>% 
  ggplot(aes(x = position, y = UC, color = group)) +
  geom_line() +
  theme_minimal() +
  labs(title = "PolyPyrimidine Probability", y = "Probability") +
  scale_color_manual(values = c("#fdb863", "#5e3c99", "#e66101", "#b2abd2")) +
  ylim(0,1)

neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_neuro_C2.fasta")
micro_neuro <- calc_pfm_long_counts(neuro_A) %>%
  mutate(group = "Neural Microexons")
nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_non_neuro_C2.fasta")
micro_nn <- calc_pfm_long_counts(nn_A) %>%
  mutate(group = "Non-Neural Microexons")
avg_nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_non_neuro_C2.fasta")
mid_nn <- calc_pfm_long_counts(avg_nn_A) %>%
  mutate(group = "Non-Neural Midexons")
avg_neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_neuro_C2.fasta")
mid_neuro <- calc_pfm_long_counts(avg_neuro_A) %>%
  mutate(group = "Neural Midexons")

polypy <- rbind(mid_neuro, mid_nn, micro_neuro, micro_nn)
compare_groups_corrected <- function(df, g1, g2) {
  df %>%
    filter(group %in% c(g1, g2)) %>%
    group_by(position) %>%
    summarise({
      table_counts <- cur_data() %>%
        group_by(group, pyrimidine) %>%
        summarise(total = sum(count), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = pyrimidine, values_from = total, values_fill = 0)
      
      mat <- as.matrix(table_counts[, c("TRUE", "FALSE")])
      rownames(mat) <- table_counts$group
      test <- fisher.test(mat)
      
      tibble(
        p_value = test$p.value,
        odds_ratio = test$estimate
      )
    }, .groups = "drop")
}

p_nmicro_navg <- compare_groups_corrected(polypy, "Neural Microexons", "Neural Midexons")
p_nmicro_nonm <- compare_groups_corrected(polypy, "Neural Microexons", "Non-Neural Microexons")
p_nmicro_nonavg <- compare_groups_corrected(polypy, "Neural Microexons", "Non-Neural Midexons")

# Merge results
p_values_df <- purrr::reduce(list(
  p_nmicro_navg, #p_value.x
  p_nmicro_nonm, #p_value.y
  p_nmicro_nonavg #p_value
), dplyr::full_join, by = "position")

# Adjust p-values for multiple comparisons 
p_values_df <- p_values_df %>%
  mutate(across(starts_with("p_value"), ~ p.adjust(.x, method = "fdr")))
p_values_df_fixed <- p_values_df %>%
  rename_with(~ paste0(., ".z"), .cols = c("p_value", "odds_ratio"))
Longer <- p_values_df_fixed %>%
  filter(position > 649,
         position < 703) %>% 
  pivot_longer(
    cols = c(odds_ratio.x, odds_ratio.y, odds_ratio.z,
             p_value.x, p_value.y, p_value.z),
    names_to = c(".value", "comparison"),
    names_pattern = "(.*)\\.(.*)"
  ) %>%
  mutate(
    comparison = recode(comparison,
                        "x" = "Micro vs Avg",
                        "y" = "Micro vs NonM",
                        "z" = "Micro vs NonAvg"),
    log2OR = log2(odds_ratio),
    signif = p_value < 0.05
  )
ggplot(Longer, aes(x = position, y = log2OR, color = signif)) +
  geom_line(aes(group = 1), linewidth = 1) +           # trend line
  geom_point(size = 2) +                          # points for each position
  geom_hline(yintercept = 0, linetype = "dashed") +  # baseline = no enrichment
  scale_color_manual(values = c("grey70", "orange")) +  # non-significant vs significant
  facet_wrap(~comparison, ncol = 1) +            # one panel per comparison
  theme_minimal(base_size = 14) +
  labs(
    x = "Nucleotide position",
    y = "log2(Odds Ratio) (Pyrimidine Enrichment)",
    color = "FDR < 0.05",
    title = "Position-wise Pyrimidine Enrichment in Microexons"
  )

################################# C15'ss #######################################
neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_neuro_C1.fasta")
# Function to calculate long-format probability matrix
calc_pfm_long <- function(seqs) {
  pfm_counts <- consensusMatrix(seqs)
  pfm_counts <- pfm_counts[c("A", "C", "G", "T"), , drop = FALSE]
  n_seqs <- length(seqs)
  df <- as.data.frame(pfm_counts)
  df$nucleotide <- rownames(df)
  rownames(df) <- NULL
  df_long <- df %>%
    pivot_longer(cols = -nucleotide,
                 names_to = "position",
                 values_to = "count"
    ) %>%
    mutate(position = as.numeric(sub("V", "", position)),
           total = n_seqs,
           probability = count / n_seqs
    )
  return(df_long)
}

micro_neuro <- calc_pfm_long(neuro_A) %>%
  mutate(group = "Neural Microexons")
summary_neuro <- calc_pfm_long(neuro_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_neuro$group <- "Neural Microexons"

nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_non_neuro_C1.fasta")
micro_nn <- calc_pfm_long(nn_A) %>%
  mutate(group = "Non-Neural Microexons")
summary_nn <- calc_pfm_long(nn_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_nn$group <- "Non-neural Microexon"

avg_nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_non_neuro_C1.fasta")
mid_nn <- calc_pfm_long(avg_nn_A) %>%
  mutate(group = "Non-Neural Midexons")
summary_avg_nn <- calc_pfm_long(avg_nn_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_avg_nn$group <- "Non-neural Midexon"

avg_neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_neuro_C1.fasta")
mid_neuro <- calc_pfm_long(avg_neuro_A) %>%
  mutate(group = "Neural Midexons")
summary_avg_neuro <- calc_pfm_long(avg_neuro_A) %>%
  group_by(nucleotide, position) %>%
  summarise(
    mean_prob = mean(probability),
    median_prob = median(probability),
    .groups = "drop"
  )
summary_avg_neuro$group <- "Neural Midexon"

mega_summary <- rbind(summary_neuro, summary_nn, summary_avg_neuro, summary_avg_nn)
mega_summary %>%
  filter(position > 200,
         position < 252) %>% 
  ggplot(aes(x = position, y = mean_prob, color = nucleotide)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Bootstrapped PFM (n = 100 each)", y = "Probability") +
  facet_grid(rows = vars(group))

## U C combo
avg_neu_polyU <- summary_avg_neuro %>%
  filter(nucleotide %in% c("T"))
avg_neu_polyU <- avg_neu_polyU[ ,c("position", "mean_prob")]
colnames(avg_neu_polyU) <- c("position", "mean_prob_T")
avg_neu_polyC <- summary_avg_neuro %>%
  filter(nucleotide %in% c("C"))
avg_neu_polyC <- avg_neu_polyC[ ,c("position", "mean_prob")]
colnames(avg_neu_polyC) <- c("position", "mean_prob_C")
avg_neu_UC <- merge(avg_neu_polyU, avg_neu_polyC, by = "position")
avg_neu_UC$UC <- avg_neu_UC$mean_prob_T + avg_neu_UC$mean_prob_C
avg_neu_UC$group <- "Neural Midexons"

avg_nn_polyU <- summary_avg_nn %>%
  filter(nucleotide %in% c("T"))
avg_nn_polyU <- avg_nn_polyU[ ,c("position", "mean_prob")]
colnames(avg_nn_polyU) <- c("position", "mean_prob_T")
avg_nn_polyC <- summary_avg_nn %>%
  filter(nucleotide %in% c("C"))
avg_nn_polyC <- avg_nn_polyC[ ,c("position", "mean_prob")]
colnames(avg_nn_polyC) <- c("position", "mean_prob_C")
avg_nn_UC <- merge(avg_nn_polyU, avg_nn_polyC, by = "position")
avg_nn_UC$UC <- avg_nn_UC$mean_prob_T + avg_nn_UC$mean_prob_C
avg_nn_UC$group <- "Non-Neural Midexons"

mic_neu_polyU <- summary_neuro %>%
  filter(nucleotide %in% c("T"))
mic_neu_polyU <- mic_neu_polyU[ ,c("position", "mean_prob")]
colnames(mic_neu_polyU) <- c("position", "mean_prob_T")
mic_neu_polyC <- summary_neuro %>%
  filter(nucleotide %in% c("C"))
mic_neu_polyC <- mic_neu_polyC[ ,c("position", "mean_prob")]
colnames(mic_neu_polyC) <- c("position", "mean_prob_C")
mic_neu_UC <- merge(mic_neu_polyU, mic_neu_polyC, by = "position")
mic_neu_UC$UC <- mic_neu_UC$mean_prob_T + mic_neu_UC$mean_prob_C
mic_neu_UC$group <- "Neural Microexons"

mic_nn_polyU <- summary_nn %>%
  filter(nucleotide %in% c("T"))
mic_nn_polyU <- mic_nn_polyU[ ,c("position", "mean_prob")]
colnames(mic_nn_polyU) <- c("position", "mean_prob_T")
mic_nn_polyC <- summary_nn %>%
  filter(nucleotide %in% c("C"))
mic_nn_polyC <- mic_nn_polyC[ ,c("position", "mean_prob")]
colnames(mic_nn_polyC) <- c("position", "mean_prob_C")
mic_nn_UC <- merge(mic_nn_polyU, mic_nn_polyC, by = "position")
mic_nn_UC$UC <- mic_nn_UC$mean_prob_T + mic_nn_UC$mean_prob_C
mic_nn_UC$group <- "Non-Neural Microexon"

polypy <- rbind(mic_neu_UC, mic_nn_UC, avg_neu_UC, avg_nn_UC)
polypy %>%
  filter(position > 200,
         position < 252) %>% 
  ggplot(aes(x = position, y = UC, color = group)) +
  geom_line() +
  theme_minimal() +
  labs(title = "PolyPyrimidine Probability", y = "Probability") +
  scale_color_manual(values = c("#fdb863", "#5e3c99", "#e66101", "#b2abd2")) +
  ylim(0,1)

neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_neuro_C1.fasta")
micro_neuro <- calc_pfm_long_counts(neuro_A) %>%
  mutate(group = "Neural Microexons")
nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/micro_image/all_non_neuro_C1.fasta")
micro_nn <- calc_pfm_long_counts(nn_A) %>%
  mutate(group = "Non-Neural Microexons")
avg_nn_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_non_neuro_C1.fasta")
mid_nn <- calc_pfm_long_counts(avg_nn_A) %>%
  mutate(group = "Non-Neural Midexons")
avg_neuro_A <- readDNAStringSet("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/avg_image/all_neuro_C1.fasta")
mid_neuro <- calc_pfm_long_counts(avg_neuro_A) %>%
  mutate(group = "Neural Midexons")

polypy <- rbind(mid_neuro, mid_nn, micro_neuro, micro_nn)
compare_groups_corrected <- function(df, g1, g2) {
  df %>%
    filter(group %in% c(g1, g2)) %>%
    group_by(position) %>%
    summarise({
      table_counts <- cur_data() %>%
        group_by(group, pyrimidine) %>%
        summarise(total = sum(count), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = pyrimidine, values_from = total, values_fill = 0)
      
      mat <- as.matrix(table_counts[, c("TRUE", "FALSE")])
      rownames(mat) <- table_counts$group
      test <- fisher.test(mat)
      
      tibble(
        p_value = test$p.value,
        odds_ratio = test$estimate
      )
    }, .groups = "drop")
}

p_nmicro_navg <- compare_groups_corrected(polypy, "Neural Microexons", "Neural Midexons")
p_nmicro_nonm <- compare_groups_corrected(polypy, "Neural Microexons", "Non-Neural Microexons")
p_nmicro_nonavg <- compare_groups_corrected(polypy, "Neural Microexons", "Non-Neural Midexons")

# Merge results
p_values_df <- purrr::reduce(list(
  p_nmicro_navg, #p_value.x
  p_nmicro_nonm, #p_value.y
  p_nmicro_nonavg #p_value
), dplyr::full_join, by = "position")

# Adjust p-values for multiple comparisons 
p_values_df <- p_values_df %>%
  mutate(across(starts_with("p_value"), ~ p.adjust(.x, method = "fdr")))
p_values_df_fixed <- p_values_df %>%
  rename_with(~ paste0(., ".z"), .cols = c("p_value", "odds_ratio"))
Longer <- p_values_df_fixed %>%
  filter(position > 200,
         position < 252) %>% 
  pivot_longer(
    cols = c(odds_ratio.x, odds_ratio.y, odds_ratio.z,
             p_value.x, p_value.y, p_value.z),
    names_to = c(".value", "comparison"),
    names_pattern = "(.*)\\.(.*)"
  ) %>%
  mutate(
    comparison = recode(comparison,
                        "x" = "Micro vs Avg",
                        "y" = "Micro vs NonM",
                        "z" = "Micro vs NonAvg"),
    log2OR = log2(odds_ratio),
    signif = p_value < 0.05
  )
ggplot(Longer, aes(x = position, y = log2OR, color = signif)) +
  geom_line(aes(group = 1), linewidth = 1) +           # trend line
  geom_point(size = 2) +                          # points for each position
  geom_hline(yintercept = 0, linetype = "dashed") +  # baseline = no enrichment
  scale_color_manual(values = c("grey70", "orange")) +  # non-significant vs significant
  facet_wrap(~comparison, ncol = 1) +            # one panel per comparison
  theme_minimal(base_size = 14) +
  labs(
    x = "Nucleotide position",
    y = "log2(Odds Ratio) (Pyrimidine Enrichment)",
    color = "FDR < 0.05",
    title = "Position-wise Pyrimidine Enrichment in Microexons"
  )
