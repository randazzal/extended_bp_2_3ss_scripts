library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(zoo)
library(ggrepel)

setwd("/data2/lackey_lab/randazza/microexons/SHAPE/")
agap1 <- read.delim("combined_5NIA/human_v_chicken_images/no_smoothing_correlation_analysis/agap_regions.txt", sep = "\t")

sliding_seq_identity <- function(df, window_size = 15) {
  # Check that required columns exist
  required_cols <- c("Seq1", "Seq2")
  if (!all(required_cols %in% names(df))) {
    stop("Dataframe must contain 'Seq1' and 'Seq2' columns.")
  }
  
  n <- nrow(df)
  identity_scores <- rep(NA, n)
  
  # Slide a window of 'window_size' across rows
  for (i in 1:(n - window_size + 1)) {
    window1 <- df$Seq1[i:(i + window_size - 1)]
    window2 <- df$Seq2[i:(i + window_size - 1)]
    
    matches <- sum(window1 == window2, na.rm = TRUE)
    identity_scores[i + floor(window_size / 2)] <- matches / window_size
  }
  
  # Add result column
  df$seq_identity <- identity_scores
  
  return(df)
}
result <- sliding_seq_identity(agap1, window_size = 15)

ggplot(result, aes(x = Nucleotide, y = seq_identity)) +
  geom_line(color = "grey60", linewidth = 0.8) +
  geom_point(
    data = subset(result, bin %in% c("anti_atandup", "cor_atandup")),
    aes(color = bin),
    size = 2.5
  ) +
  scale_color_manual(
    values = c("anti_atandup" = "#E64B35", "cor_atandup" = "#4DBBD5"),
    name = "Highlighted bins"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Nucleotide position",
    y = "Sequence identity",
    title = "Sliding 23-nt sequence identity with highlighted bins"
  )
highlight_regions <- result |>
  dplyr::filter(bin %in% c("cor_atandup", "cor_below")) |>
  dplyr::group_by(bin, grp = cumsum(c(1, diff(Nucleotide) != 1))) |>
  dplyr::summarize(xmin = min(Nucleotide), xmax = max(Nucleotide), .groups = "drop")

ggplot(result, aes(x = Nucleotide, y = seq_identity)) +
  geom_rect(
    data = highlight_regions,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = bin),
    alpha = 0.4, inherit.aes = FALSE
  ) +
  geom_line(color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("cor_atandup" = "#b2df8a", "cor_below" = "#a6cee3")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Nucleotide position",
    y = "Sequence identity",
    title = "Sliding sequence identity with shaded bins"
  ) +
  geom_vline(xintercept = c(266, 277)) ##first and last position of microexon (see below)

#agap1 266:277; apbb2 283:288; asap2 265:273; clec16a 264:269; cpeb4 287:310; dctn1 289:304; dock7 281:295; 
#emc1 279:288; fry 273:281; itsn1 279:294; mef2a 275:300; ptprk 280:292; robo1 278:287
