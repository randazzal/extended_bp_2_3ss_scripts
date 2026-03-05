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
#################################### determining length cutoffs for regions of interest
# --- Function: rolling correlations, and run stats ---
shuffle_run_distribution <- function(df1, df2, 
                                     n_iter = 1000,
                                     corr_window = 15,
                                     corr_method = "spearman",
                                     seed = 123) {
  set.seed(seed)
  
  half_corr <- floor(corr_window / 2)
  min_pairs <- ceiling(corr_window * 0.5)  # require at least half valid points
  
  # Helper: rolling correlation for one sequence pair
  win_cor <- function(x, y, method = corr_method) {
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < min_pairs) return(NA)
    cor(x[ok], y[ok], method = method)
  }
  
  # Store all run lengths
  all_runs <- numeric()
  
  for (i in seq_len(n_iter)) {
    if (i %% 50 == 0) cat("Iteration", i, "of", n_iter, "\n")
    
    # Shuffle reactivities
    shuffled_df1 <- df1 %>% mutate(React1 = sample(React1))
    shuffled_df2 <- df2 %>% mutate(React2 = sample(React2))
    
    # cap at [0,4]
    shuffled_df1 <- shuffled_df1 %>%
      mutate(React1 = pmin(pmax(React1, 0), 4))
    shuffled_df2 <- shuffled_df2 %>%
      mutate(React2 = pmin(pmax(React2, 0), 4))
    
    # Merge and compute rolling correlation
    merged <- inner_join(shuffled_df1, shuffled_df2, by = "Nucleotide")
    mat <- cbind(merged$React1, merged$React2)
    
    rho <- rollapply(
      data = mat,
      width = corr_window,
      by = 1,
      FUN = function(m) win_cor(m[,1], m[,2], method = corr_method),
      by.column = FALSE,
      align = "center",
      partial = FALSE
    )
    
    # Pad edges
    res <- data.frame(
      position = seq_along(merged$React1),
      rho = c(rep(NA, half_corr), rho, rep(NA, half_corr))
    )
    
    # Find runs of positive correlation
    pos_runs <- rle(res$rho < -0.2)
    run_lengths <- pos_runs$lengths[pos_runs$values == TRUE]
    
    all_runs <- c(all_runs, run_lengths)
  }
  
  # Summarize results
  run_distribution <- as.data.frame(table(all_runs))
  run_distribution <- run_distribution %>%
    mutate(all_runs = as.numeric(as.character(all_runs))) %>%
    arrange(all_runs)
  
  list(
    run_lengths = all_runs,
    distribution = run_distribution
  )
}

agap1 <- read.delim("combined_5NIA/human_v_chicken_images/trimmed_AGAP1_aligned_reactivities.txt")
agap1$Nucleotide <- 1:498
agap1_hsa <- agap1[, c("React1", "Nucleotide")]
agap1_gga <- agap1[, c("React2", "Nucleotide")]

result_runs <- shuffle_run_distribution(
  agap1_hsa,
  agap1_gga,
  n_iter = 1000,
  corr_window = 15,
  corr_method = "spearman",
  seed = 123
)

head(result_runs$distribution)
hist(result_runs$run_lengths, breaks = 30,
     main = "Distribution of All >0.2 rho Run Lengths",
     xlab = "Run length (consecutive rho > 0.2)")
sum(result_runs$distribution$Freq[result_runs$distribution$all_runs >= 19]) /
  sum(result_runs$distribution$Freq) * 100

########################### unshuffled rho runs
positive_rho_run_distribution <- function(df1, df2, 
                                          corr_window = 15,
                                          corr_method = "spearman",
                                          threshold = 0.2,
                                          run_cutoff = 19) {
  
  # Helper: rolling correlation
  min_pairs <- ceiling(corr_window * 0.5)
  win_cor <- function(x, y, method = corr_method) {
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < min_pairs) return(NA)
    cor(x[ok], y[ok], method = method)
  }
  
  # Preprocess: replace NA with 0, cap at [0,4], compute rolling averages
  df1 <- df1 %>%
    mutate(React1 = pmin(pmax(React1, 0), 4))
  df2 <- df2 %>%
    mutate(React2 = pmin(pmax(React2, 0), 4))
  
  # Merge by position
  merged <- inner_join(df1, df2, by = "Nucleotide")
  mat <- cbind(merged$React1, merged$React2)
  
  # Rolling correlation
  half_corr <- floor(corr_window / 2)
  rho <- rollapply(
    data = mat,
    width = corr_window,
    by = 1,
    FUN = function(m) win_cor(m[,1], m[,2], method = corr_method),
    by.column = FALSE,
    align = "center",
    partial = FALSE
  )
  
  # Pad edges
  res <- data.frame(
    Nucleotide = merged$Nucleotide,
    rho = c(rep(NA, half_corr), rho, rep(NA, half_corr))
  )
  
  # Categorize correlations
  res$category <- case_when(
    res$rho >= threshold ~ "correlated",
    res$rho <= -threshold ~ "anti_correlated",
    !is.na(res$rho) ~ "neutral",
    TRUE ~ NA_character_
  )
  
  # Run length encoding
  runs <- rle(res$category)
  
  # Initialize bin column
  res$bin <- NA_character_
  
  # Assign bin labels for correlated and anti-correlated runs
  idx <- 1
  for (i in seq_along(runs$lengths)) {
    len <- runs$lengths[i]
    val <- runs$values[i]
    
    if (is.na(val)) {
      idx <- idx + len
      next
    }
    
    if (val == "correlated") {
      bin_label <- if (len < run_cutoff) "cor_below" else "cor_atandup"
      res$bin[idx:(idx + len - 1)] <- bin_label
    } else if (val == "anti_correlated") {
      bin_label <- if (len < run_cutoff) "anti_below" else "anti_atandup"
      res$bin[idx:(idx + len - 1)] <- bin_label
    }
    
    idx <- idx + len
  }
  
  # Run lengths by category
  correlated_runs <- runs$lengths[runs$values == "correlated"]
  anti_runs <- runs$lengths[runs$values == "anti_correlated"]
  neutral_runs <- runs$lengths[runs$values == "neutral"]
  
  # Return everything
  list(
    correlated = correlated_runs,
    neutral = neutral_runs,
    anti_correlated = anti_runs,
    rho_table = table(res$category, useNA = "ifany"),
    rho_df = res
  )
}
result_runs <- positive_rho_run_distribution(agap1_hsa, agap1_gga,
                                             corr_window = 15,
                                             corr_method = "spearman",
                                             threshold = 0.2,
                                             run_cutoff = 19)

# View counts for each category
result_runs$rho_table
# Distribution of correlated run lengths
table(result_runs$correlated)
table(result_runs$anti_correlated)

agap_start <- merge(agap1, result_runs$rho_df, by = "Nucleotide")
write.table(agap_start, "combined_5NIA/human_v_chicken_images/no_smoothing_correlation_analysis/agap_regions.txt", sep = "\t")
