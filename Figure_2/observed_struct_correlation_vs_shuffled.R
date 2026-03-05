library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(zoo)
library(tidyverse)
library(ggbeeswarm)
library(tibble)

setwd("/data2/lackey_lab/randazza/microexons/SHAPE/")
# --- Function to calculate rolling correlation without plotting ---
rolling_correlation <- function(df1, df2, method = "spearman", return_merged = TRUE) {
  # Replace NA with 0 and cap values
  df1 <- df1 %>% mutate(React1 = pmin(pmax(React1, 0), 4))
  df2 <- df2 %>% mutate(React2 = pmin(pmax(React2, 0), 4))
  # Merge by position
  merged <- inner_join(df1, df2, by = "Nucleotide")
  # Compute correlation
  corr_val <- cor(merged$React1, merged$React2, use = "complete.obs", method = method)
  # Return rho
  if (return_merged) {
    return(list(rho = corr_val, merged = merged))
  } else {
    return(corr_val)
  }
}
# --- Shuffle test ---
shuffle_correlations <- function(df1, df2, method = "spearman",
                                 n_iter = 1000, seed = 123, save_merged_iter = 7) {
  set.seed(seed)
  r_values <- numeric(n_iter)
  merged_output <- NULL
  
  for (i in 1:n_iter) {
    shuffled_df1 <- df1 %>% mutate(React1 = sample(React1))
    shuffled_df2 <- df2 %>% mutate(React2 = sample(React2))
    
    if (i == save_merged_iter) {
      # Save merged dataframe for this iteration
      result <- rolling_correlation(shuffled_df1, shuffled_df2, method, return_merged = TRUE)
      r_values[i] <- result$rho
      merged_output <- result$merged
    } else {
      r_values[i] <- rolling_correlation(shuffled_df1, shuffled_df2, method)
    }
  }
  
  list(
    r_df = tibble(iteration = 1:n_iter, rho = as.numeric(r_values)),
    merged_example = merged_output
  )
}
agap1 <- read.delim("combined_5NIA/human_v_chicken_images/trimmed_AGAP1_aligned_reactivities.txt")
agap1$Nucleotide <- 1:498
agap1_hsa <- agap1[ ,c("React1", "Nucleotide")]
agap1_gga <- agap1[ ,c("React2", "Nucleotide")]
result <- shuffle_correlations(agap1_hsa, agap1_gga, n_iter = 1000, save_merged_iter = 98)
r_results <- result$r_df 
r_results$group <- "Shuffled AGAP1"
summary(r_results$rho)
ggplot(r_results, aes(x = rho)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(
    title = "Null distribution of Correlations (shuffled reactivities)",
    x = "Spearman Correlation",
    y = "Count"
  ) +
  geom_vline(xintercept = 0.266)
df <- data.frame(
  gene = c("ROBO1", "PTPRK", "FRY", "CLEC16A", "APBB2", "ITSN1", "AGAP1", "EMC1", "DOCK7", "DCTN1", "ASAP2", "MEF2A", "CPEB4"),
  rho = as.numeric(c("0.47", "0.272", "0.269", "0.262", "0.327", "0.278", "0.392", "0.085", "0.451", "-0.048", "0.225", "0.285", "0.193"))
)
df$group <- "Observed"
combo <- rbind(r_results[ ,c("rho", "group")], df[ ,c("rho", "group")])
combo$rho <- as.numeric(combo$rho)
combo %>%
  ggplot(aes(x = group, y = rho, fill = group)) +
  geom_boxplot() +
  geom_quasirandom(width = 0.15, size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Observed vs Shuffled Spearman Distribution",
    x = "",
    y = "rho"
  ) +
  scale_y_log10() +
  ggeasy::easy_center_title()

emc1 <- read.delim("combined_5NIA/human_v_chicken_images/trimmed_EMC1_aligned_reactivities.txt")
emc1$Nucleotide <- 1:531
emc1_hsa <- emc1[ ,c("React1", "Nucleotide")]
emc1_gga <- emc1[ ,c("React2", "Nucleotide")]
result <- shuffle_correlations(emc1_hsa, emc1_gga, n_iter = 1000, save_merged_iter = 98)
r_results_emc <- result$r_df 
r_results_emc$group <- "Shuffled EMC1"
summary(r_results_emc$rho)
ggplot(r_results_emc, aes(x = rho)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(
    title = "Null distribution of Correlations (shuffled reactivities)",
    x = "Spearman Correlation",
    y = "Count"
  ) +
  geom_vline(xintercept = 0.266)

combo <- rbind(r_results[ ,c("rho", "group")], r_results_emc[ ,c("rho", "group")], df[ ,c("rho", "group")])
combo$rho <- as.numeric(combo$rho)
combo %>%
  ggplot(aes(x = group, y = rho, fill = group)) +
  geom_boxplot() +
  geom_quasirandom(width = 0.15, size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Observed vs Shuffled Spearman Distribution",
    x = "",
    y = "rho"
  ) +
  scale_y_log10() +
  ggeasy::easy_center_title()
combo %>%
  ggplot(aes(x = group, y = rho, fill = group)) +
  geom_boxplot() +
  geom_quasirandom(width = 0.15, size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Observed vs Shuffled Spearman Distribution",
    x = "",
    y = "rho"
  ) +
  ggeasy::easy_center_title() +
  scale_fill_manual(values = c("Observed" = "#b3e2cd", "Shuffled AGAP1" = "#cbd5e8", "Shuffled EMC1" = "#f4cae4"))

wilcox.test(r_results$rho, r_results_emc$rho) #0.9762
wilcox.test(r_results$rho, df$rho) #8.074e-08
wilcox.test(r_results_emc$rho, df$rho) #7.905e-08